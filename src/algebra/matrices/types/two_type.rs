

use crate::{algebra::{matrices::{
    operations::umatch::differential::DifferentialUmatch,
    query::{MatrixAlgebra, MatrixOracle},
    types::transpose::OrderAntiTranspose,
}, rings::traits::SemiringOperations, vectors::entries::KeyValGet}, utilities::{iterators::general::TwoTypeIterator, order::TwoTypeOrderOperator}};


use std::hash::Hash;
use std::vec::IntoIter;
use std::fmt::Debug;

pub enum TwoTypeMatrix< Matrix1, Matrix2 >{
    Version1( Matrix1 ),
    Version2( Matrix2 ),
}



// -----------------------------------------
// MATRIX ORACLE
// -----------------------------------------


impl <
    Matrix1,
    Matrix2,
    RowIndex,
    ColumnIndex,
    RowEntry,
    ColumnEntry,
    Coefficient,
>

    MatrixOracle for

    TwoTypeMatrix
        < Matrix1, Matrix2 >

    where
        Matrix1:            MatrixOracle<
                                RowIndex = RowIndex,
                                ColumnIndex = ColumnIndex,
                                RowEntry = RowEntry,
                                ColumnEntry = ColumnEntry,
                                Coefficient = Coefficient,
                            >,
        Matrix2:            MatrixOracle<
                                RowIndex = RowIndex,
                                ColumnIndex = ColumnIndex,
                                RowEntry = RowEntry,
                                ColumnEntry = ColumnEntry,
                                Coefficient = Coefficient,
                            >,

        // these are type constraints imposed on anything that implements the `MatrixOracle` trait
        Coefficient:            Clone + Debug + PartialEq,    // The type of coefficient stored in each entry of the matrix    
        RowIndex   :            Clone + Debug + Eq,    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
        ColumnIndex:            Clone + Debug + Eq,    // The type of column indices    
        RowEntry   :            Clone + Debug + PartialEq + KeyValGet <  Key = ColumnIndex,  Val = Coefficient   >,  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
        ColumnEntry:            Clone + Debug + PartialEq + KeyValGet <  Key = RowIndex,     Val = Coefficient   >,  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`

{
    type Coefficient            =   Coefficient;    // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   RowIndex;    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   ColumnIndex;    // The type of column indices    
    type RowEntry               =   RowEntry;    // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   ColumnEntry;  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =   TwoTypeIterator< Matrix1::Row, Matrix2::Row >;  // What you get when you ask for a row
    type RowReverse             =   TwoTypeIterator< Matrix1::RowReverse, Matrix2::RowReverse >;  // What you get when you ask for a row with the order of entries reversed
    type Column                 =   TwoTypeIterator< Matrix1::Column, Matrix2::Column >;  // What you get when you ask for a column
    type ColumnReverse          =   TwoTypeIterator< Matrix1::ColumnReverse, Matrix2::ColumnReverse >;  // What you get when you ask for a column with the order of entries reversed



    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeIterator::Version1(matrix.row(index)),
            TwoTypeMatrix::Version2(matrix) => TwoTypeIterator::Version2(matrix.row(index)),
        }
    }
    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeIterator::Version1(matrix.row_reverse(index)),
            TwoTypeMatrix::Version2(matrix) => TwoTypeIterator::Version2(matrix.row_reverse(index)),
        }
    }
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeIterator::Version1(matrix.column(index)),
            TwoTypeMatrix::Version2(matrix) => TwoTypeIterator::Version2(matrix.column(index)),
        }
    }
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeIterator::Version1(matrix.column_reverse(index)),
            TwoTypeMatrix::Version2(matrix) => TwoTypeIterator::Version2(matrix.column_reverse(index)),
        }
    }
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        match self {
            TwoTypeMatrix::Version1(matrix) => matrix.has_row_for_index(index),
            TwoTypeMatrix::Version2(matrix) => matrix.has_row_for_index(index),
        }
    }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        match self {
            TwoTypeMatrix::Version1(matrix) => matrix.has_column_for_index(index),
            TwoTypeMatrix::Version2(matrix) => matrix.has_column_for_index(index),
        }
    }


    /// Returns `Some(x)` if there is a structural nonzero at `(row,column)`. 
    /// 
    ///  Returns `None` otherwise.
    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        match self {
            TwoTypeMatrix::Version1(matrix) => matrix.structural_nonzero_entry(row, column),
            TwoTypeMatrix::Version2(matrix) => matrix.structural_nonzero_entry(row, column),
        }
    }
  
}






// -----------------------------------------
// MATRIX ALGEBRA
// -----------------------------------------




impl <
    Matrix1,
    Matrix2,
    RowIndex,
    ColumnIndex,
    RowEntry,
    ColumnEntry,
    Coefficient,
    RingOperator,
>

    MatrixAlgebra for

    TwoTypeMatrix
        < Matrix1, Matrix2 >

    where
        Matrix1:            MatrixAlgebra<
                                RowIndex = RowIndex,
                                ColumnIndex = ColumnIndex,
                                RowEntry = RowEntry,
                                ColumnEntry = ColumnEntry,
                                Coefficient = Coefficient,
                                RingOperator = RingOperator,
                            >,
        Matrix2:            MatrixAlgebra<
                                RowIndex = RowIndex,
                                ColumnIndex = ColumnIndex,
                                RowEntry = RowEntry,
                                ColumnEntry = ColumnEntry,
                                Coefficient = Coefficient,
                                RingOperator = RingOperator,
                            >,

        // these are type constraints imposed on anything that implements the `MatrixOracle` trait
        Coefficient:            Clone + Debug + PartialEq,    // The type of coefficient stored in each entry of the matrix    
        RowIndex   :            Clone + Debug + Eq,    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
        ColumnIndex:            Clone + Debug + Eq,    // The type of column indices    
        RowEntry   :            Clone + Debug + PartialEq + KeyValGet <  Key = ColumnIndex,  Val = Coefficient   >,  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
        ColumnEntry:            Clone + Debug + PartialEq + KeyValGet <  Key = RowIndex,     Val = Coefficient   >,  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`


        RingOperator:                   Clone + SemiringOperations< Element = Self::Coefficient >,
        // OrderOperatorForRowEntries:     Clone + JudgeOrder< Self::RowEntry >,
        // OrderOperatorForRowIndices:     Clone + JudgeOrder< Self::RowIndex >,
        // OrderOperatorForColumnEntries:  Clone + JudgeOrder< Self::ColumnEntry >,
        // OrderOperatorForColumnIndices:  Clone + JudgeOrder< Self::ColumnIndex >,


{
    type RingOperator = RingOperator; // The type of ring operator used to perform arithmetic operations on the coefficients of the matrix entries  

    type OrderOperatorForRowEntries = TwoTypeOrderOperator<
        Matrix1::OrderOperatorForRowEntries,
        Matrix2::OrderOperatorForRowEntries
        >;

    type OrderOperatorForRowIndices = TwoTypeOrderOperator<
        Matrix1::OrderOperatorForRowIndices,
        Matrix2::OrderOperatorForRowIndices
    >;

    type OrderOperatorForColumnEntries = TwoTypeOrderOperator<
        Matrix1::OrderOperatorForColumnEntries,
        Matrix2::OrderOperatorForColumnEntries
    >;

    type OrderOperatorForColumnIndices = TwoTypeOrderOperator<
        Matrix1::OrderOperatorForColumnIndices,
        Matrix2::OrderOperatorForColumnIndices
    >;

    fn ring_operator( &self ) -> Self::RingOperator {
        match self {
            TwoTypeMatrix::Version1(matrix) => matrix.ring_operator(),
            TwoTypeMatrix::Version2(matrix) => matrix.ring_operator(),
        }
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeOrderOperator::Version1(matrix.order_operator_for_row_entries()),
            TwoTypeMatrix::Version2(matrix) => TwoTypeOrderOperator::Version2(matrix.order_operator_for_row_entries()),
        }
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeOrderOperator::Version1(matrix.order_operator_for_row_indices()),
            TwoTypeMatrix::Version2(matrix) => TwoTypeOrderOperator::Version2(matrix.order_operator_for_row_indices()),
        }
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeOrderOperator::Version1(matrix.order_operator_for_column_entries()),
            TwoTypeMatrix::Version2(matrix) => TwoTypeOrderOperator::Version2(matrix.order_operator_for_column_entries()),
        }
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        match self {
            TwoTypeMatrix::Version1(matrix) => TwoTypeOrderOperator::Version1(matrix.order_operator_for_column_indices()),
            TwoTypeMatrix::Version2(matrix) => TwoTypeOrderOperator::Version2(matrix.order_operator_for_column_indices()),
        }
    }
}        




























#[cfg(test)]
mod test {
    use itertools::Itertools;

    use super::*;
    use crate::algebra::{matrices::{debug::{matrix_oracle_is_internally_consistent, matrix_order_operators_are_internally_consistent, product_is_identity_matrix}, types::{packet::MatrixAlgebraPacket, product::ProductMatrix, vec_of_vec::sorted::VecOfVec}}, rings::types::field_prime_order::PrimeOrderField};












    /// Checks that the `TwoTypeMatrix` is internally consistent.
    #[test]
    fn comprehensive_test() {   

        use crate::algebra::matrices::operations::umatch::row_major::Umatch;
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::query::MatrixOracle;               

        let num_indices_row             =   10;
        let num_indices_col             =   20;
        let approximate_density           =   0.2;
        let modulus                     =   17;
        let allow_nonstructural_zero     =   true;

        let ring_operator           =   PrimeOrderField::new( modulus );
        let matrix_data       =   VecOfVec::random_mod_p_with_density( num_indices_row, num_indices_col, approximate_density, modulus, allow_nonstructural_zero );
        let matrix = MatrixAlgebraPacket::with_default_order( & matrix_data, ring_operator );
     

        let row_indices = (0..num_indices_row).rev().collect_vec();
        let column_indices = (0..num_indices_col).collect_vec();

        let get_matrix = | version1: bool | {
            if version1 {
                TwoTypeMatrix::Version1(
                     MatrixAlgebraPacket::with_default_order(
                        matrix.clone(),
                        ring_operator.clone()
                    )
                )
            } else {
                TwoTypeMatrix::Version2( 
                    MatrixAlgebraPacket::with_default_order(
                        matrix.clone(),
                        ring_operator.clone()
                    ) 
                )
            }
        };

        let a = get_matrix( true );
        let b = get_matrix( false );

        // check matrices are internally consistent

        assert!(
            matrix_oracle_is_internally_consistent(
                & a, 
                0..num_indices_row, 
                0..num_indices_col
            )
        );

        assert!(
            matrix_oracle_is_internally_consistent(
                & b, 
                0..num_indices_row, 
                0..num_indices_col
            )
        );   

        // check order operators are correct        

        assert!(
            matrix_order_operators_are_internally_consistent(
                & a, 
                0..num_indices_row, 
                0..num_indices_col
            ).is_ok()
        );

        assert!(
            matrix_order_operators_are_internally_consistent(
                b, 
                0..num_indices_row, 
                0..num_indices_col
            ).is_ok()
        );     
              
    }







}