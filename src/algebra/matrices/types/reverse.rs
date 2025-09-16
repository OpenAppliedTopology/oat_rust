//! Lazy transpose and anti-transpose; wraps around another matrix, and swaps order and/or major vs. columns.

use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle};
use crate::algebra::matrices::operations::{MatrixOracleOperations};
use crate::utilities::order::ReverseOrder;






//  REVERSE
//  -----------------------------------------------------------------------------------------------

/// Wraps around a matrix; returns entries of each row (respectively, column) in reverse order (warning: doesn't work quite the same as for ordinary dense matrices)
/// 
/// Concretely, the antitranpose is obtained by (i) transposing, and then (ii) reversing the order of rows and reversing the order of columns.
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  For example, to execute the method `row`, an `Reverse` struct simply calls `row_reverse` on
/// the underlying matrix oracle.   **Note that this doesn't work quite the same as for ordinary dense matrices when indices happen to be integers.**
/// 
/// # Caution
/// 
/// There are three important differences between [Reverse] and the matrix returned by [antitranspose_deep](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep):
/// 
/// - [Reverse] is a lazy object that does not generate new data
/// - The set of (key,val) pairs that appear in a major (respectively, minor) view of `Reverse::new(matrix)`
///   are the *same* as the entries in a minor (respectively, major) view of `matrix`; only the sequence in which those entries appear is different.
///   By contrast, the keys in the (key,val) pairs of [matrix.antitranspose_deep()](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) are different;
///   they are obtained by subtracting the original keys from (# rows in the antitransposed matrix - 1).
/// - For this reason, [matrix.antitranspose_deep()](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) is only available for
///   very specific types of matrices; [Reverse] is available for a much broader class.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::{reverse::ReverseMatrix, transpose::Transpose};
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::utilities::iterators::general::result_iterators_are_elementwise_equal;
/// use itertools;
/// 
/// // matrix
/// let matrix =   & VecOfVec::new( vec![
///                                     vec![ (0,0), (1,1), (2,2) ],
///                                     vec![ (0,3), (1,4), (2,5) ],
///                                 ] ).ok().unwrap();
/// 
/// // its reverse
/// let rev    =   ReverseMatrix::new( matrix );
/// 
/// // check that rows are reversed correctly
/// for index in 0 .. 2 {
///     assert!( matrix.row_reverse(&index).eq(        rev.row(&index)              )           );
///     assert!( matrix.row(&index).eq(                rev.row_reverse(&index)      )           );
/// }
/// 
/// // check that columns are reversed correctly
/// for index in 0 .. 3 {
///     assert!( matrix.column(&index).eq(             rev.column_reverse(&index)   )           );
///     assert!( matrix.column_reverse(&index).eq(     rev.column(&index)           )           );
/// }  
/// 
/// // check that valid and invalid indices return the correct results
/// for index in 0 .. 5 {    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_result(&index),                  rev.row_reverse_result(&index) 
///     ));    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_reverse_result(&index),          rev.row_result(&index) 
///     ));    
/// 
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_result(&index),               rev.column_reverse_result(&index) 
///     ));    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_reverse_result(&index),       rev.column_result(&index) 
///     ));                                     
/// }
/// ```
pub struct ReverseMatrix< Matrix > { unreversed_matrix: Matrix }

impl < Matrix >

    ReverseMatrix
        < Matrix >
{
    pub fn new( unreversed_matrix: Matrix ) -> ReverseMatrix< Matrix > { ReverseMatrix { unreversed_matrix: unreversed_matrix } }
}


    //  CLONE
impl < Matrix: Clone > 

    Clone for

    ReverseMatrix< Matrix >

{ fn clone(&self) -> Self { ReverseMatrix { unreversed_matrix: self.unreversed_matrix.clone() } } }    


//  MATRIX ORACLE
//  ---------------------------------------------------------------------------


impl< Matrix > 

    MatrixOracle for 
    
    ReverseMatrix< Matrix >

    where
        Matrix:     MatrixOracle
{
    type Coefficient        =   Matrix::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   Matrix::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   Matrix::ColumnIndex;       // The type of column indices
    
    type RowEntry           =   Matrix::RowEntry;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry        =   Matrix::ColumnEntry;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                =   Matrix::RowReverse;               // What you get when you ask for a row.
    type RowReverse         =   Matrix::Row;        // What you get when you ask for a row with the order of entries reversed
    
    type Column             =   Matrix::ColumnReverse;            // What you get when you ask for a column
    type ColumnReverse      =   Matrix::Column;     // What you get when you ask for a column with the order of entries reversed 

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        self.unreversed_matrix.structural_nonzero_entry( row, column )
    }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool
        { self.unreversed_matrix.has_row_for_index( index ) }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool
        { self.unreversed_matrix.has_column_for_index( index ) }    

    fn row(                     &self,  index: &Self::RowIndex    )       -> Self::Row 
        { self.unreversed_matrix.row_reverse(index) }
    fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex >
        { self.unreversed_matrix.row_reverse_result(index) }     
    fn row_reverse(             &self,  index: &Self::RowIndex    )       -> Self::RowReverse
        { self.unreversed_matrix.row(index) }    
    fn row_reverse_result(         &   self, index: & Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex >
        { self.unreversed_matrix.row_result(index) }
    
    fn column(                  &self,  index: &Self::ColumnIndex )       -> Self::Column
        { self.unreversed_matrix.column_reverse(index) }
    fn column_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex >
        { self.unreversed_matrix.column_reverse_result(index) }    
    fn column_reverse(          &self,  index: &Self::ColumnIndex )       -> Self::ColumnReverse
        { self.unreversed_matrix.column(index) }
    fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex >
        { self.unreversed_matrix.column_result(index) }
} 





//  MATRIX ALGEBRA
//  ---------------------------------------------------------------------------


impl< Matrix > 

    MatrixAlgebra for 
    
    ReverseMatrix< Matrix >

    where
        Matrix:     MatrixAlgebra
{
    type RingOperator                                   =   Matrix::RingOperator;

    type OrderOperatorForRowEntries                     =   ReverseOrder< Matrix::OrderOperatorForRowEntries >;

    type OrderOperatorForRowIndices                     =   ReverseOrder< Matrix::OrderOperatorForRowIndices >;

    type OrderOperatorForColumnEntries                  =   ReverseOrder< Matrix::OrderOperatorForColumnEntries >;

    type OrderOperatorForColumnIndices                  =   ReverseOrder< Matrix::OrderOperatorForColumnIndices >;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.unreversed_matrix.ring_operator()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        ReverseOrder::new(
            self.unreversed_matrix.order_operator_for_row_entries()
        )
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        ReverseOrder::new(
            self.unreversed_matrix.order_operator_for_row_indices()
        )
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        ReverseOrder::new(
            self.unreversed_matrix.order_operator_for_column_entries()
        )
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        ReverseOrder::new(
            self.unreversed_matrix.order_operator_for_column_indices()
        )
    }
}




//  MATRIX ORACLE OPERATIONS
//  ---------------------------------------------------------------------------


impl< Matrix > 

    MatrixOracleOperations for 
    
    ReverseMatrix< Matrix >
{}    