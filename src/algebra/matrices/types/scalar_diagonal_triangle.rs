//! Matrix wrapper that prepends an entry to each row and postpends an entry to each column of the wrapped matrix.


use std::{iter::{Once, once}};

use crate::algebra::matrices::{operations::MatrixOracleOperations, query::{MatrixAlgebra, MatrixOracle}};
use crate::algebra::vectors::entries::KeyValNew;



/// Represents the sum of a scalar matrix and a strictly upper triangular matrix.
/// 
/// Recall that a scalar matrix is a diagonal matrix where all diagonal entries are equal.
/// A matrix `M` is strictly upper triangular if `i < j` whenever `M[i,j]` is nonzero.
/// 
/// This data structure is used to save on memory, because matrices of this form are common
/// in applied topology, and it's often the case that the strictly upper triangular part has
/// very few nonzero entries.  Thus it becomes efficient to store one single scalar for the
/// whole diagonal, and a small number of off-diagonal entries.
/// 
/// # Warning
/// 
/// This data structure is unsafe, in the sense that it does not check that the strictly-upper-trianuglar matrix
/// that the user provides is actually strictly upper triangular. It simply appends an element representing the
/// the diagonal element to the beginning of each row vector and reversed column vector (respectively, to the end of each reversed row vector and each column vector), and returns the user-provided value for
/// the diagonal element whenever the user asks for an element anywhere on the diagonal.
/// 
/// # Design notes
/// 
/// The developers originally attempted to write this struct such that the `strictly_upper_triangular_matrix`
/// field *owned* `StrictlyUpperTriangularMatrix`, not just a *reference* to `StrictlyUpperTriangularMatrix`.  However,
/// this lead to errors which we didn't fully understand.  We hypothesize that this may have to
/// do with lifetime parameters.  In particular, in the original version where `strictly_upper_triangular_matrix`
/// owned its matrix and not a reference, the compiler asked us to give an explicit type
/// annotation for the `Row` type of `StrictlyUpperTriangularMatrix`; that ran us into difficulty,
/// becasue the `Row` type involved a lifetime parameter.
pub struct SumOfScalarAndStrictlyUpperTriangularMatrices < StrictlyUpperTriangularMatrix >
    where 
        StrictlyUpperTriangularMatrix:          MatrixOracle
{
    strictly_upper_triangular_matrix:           StrictlyUpperTriangularMatrix,
    diagonal_scalar:                            StrictlyUpperTriangularMatrix::Coefficient,
}

// Implement the struct
impl < StrictlyUpperTriangularMatrix: MatrixOracle >

    SumOfScalarAndStrictlyUpperTriangularMatrices 
        < StrictlyUpperTriangularMatrix >
{
    pub fn new( strictly_upper_triangular_matrix: StrictlyUpperTriangularMatrix, diagonal_scalar: StrictlyUpperTriangularMatrix::Coefficient ) -> Self {  
        SumOfScalarAndStrictlyUpperTriangularMatrices{
                strictly_upper_triangular_matrix,
                diagonal_scalar, 
            }
    }

}



//  Matrix Oracle
//  ---------------------------------------------------------------------------


impl < StrictlyUpperTriangularMatrix, Index > 

    MatrixOracle for 

    SumOfScalarAndStrictlyUpperTriangularMatrices 
        < StrictlyUpperTriangularMatrix >

    where 
        StrictlyUpperTriangularMatrix:      MatrixOracle< 
                                                RowIndex        =   Index,          // we use the Index type parameter to ensure that row and column indices to have the same type
                                                ColumnIndex     =   Index,          // we use the Index type parameter to ensure that row and column indices to have the same type
                                                RowEntry:           KeyValNew,
                                                ColumnEntry:        KeyValNew
                                            >, 
        Index:                              Clone + Eq + std::fmt::Debug,
{
    type Coefficient            =   StrictlyUpperTriangularMatrix::Coefficient   ;

    type RowIndex               =   StrictlyUpperTriangularMatrix::RowIndex      ;

    type ColumnIndex            =   StrictlyUpperTriangularMatrix::ColumnIndex   ;

    type RowEntry               =   StrictlyUpperTriangularMatrix::RowEntry      ;

    type ColumnEntry            =   StrictlyUpperTriangularMatrix::ColumnEntry   ;

    type Row                    =   core::iter::Chain < Once < StrictlyUpperTriangularMatrix::RowEntry >, StrictlyUpperTriangularMatrix::Row >;

    type RowReverse             =   core::iter::Chain < StrictlyUpperTriangularMatrix::RowReverse, Once < StrictlyUpperTriangularMatrix::RowEntry > >;

    type Column                 =   core::iter::Chain < StrictlyUpperTriangularMatrix::Column, Once < StrictlyUpperTriangularMatrix::ColumnEntry > >;

    type ColumnReverse          =   core::iter::Chain < Once < StrictlyUpperTriangularMatrix::ColumnEntry >, StrictlyUpperTriangularMatrix::ColumnReverse >;

    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        std::iter::Iterator::chain(
            once( Self::RowEntry::new( index.clone(), self.diagonal_scalar.clone() )  ),
            self.strictly_upper_triangular_matrix.row( & index ),
        )
    }

    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        std::iter::Iterator::chain(
            self.strictly_upper_triangular_matrix.row_reverse( & index ),
            once( Self::RowEntry::new( index.clone(), self.diagonal_scalar.clone() )  ),            
        )
    }

    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
        std::iter::Iterator::chain(
            self.strictly_upper_triangular_matrix.column( & index ),
            once( Self::ColumnEntry::new( index.clone(), self.diagonal_scalar.clone() )  ),            
        )
    }

    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        std::iter::Iterator::chain(
            once( Self::ColumnEntry::new( index.clone(), self.diagonal_scalar.clone() )  ),            
            self.strictly_upper_triangular_matrix.column_reverse( & index ),            
        )
    }

    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.strictly_upper_triangular_matrix.has_row_for_index( index )
    }

    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.strictly_upper_triangular_matrix.has_column_for_index( index )
    }

    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        match row == column {
            true => { return Some( self.diagonal_scalar.clone() ) },
            false => { return self.strictly_upper_triangular_matrix.structural_nonzero_entry( row, column ) }
        }
    }
}    







//  Matrix Algebra
//  ---------------------------------------------------------------------------


impl < StrictlyUpperTriangularMatrix, Index > 

    MatrixAlgebra for 

    SumOfScalarAndStrictlyUpperTriangularMatrices 
        < StrictlyUpperTriangularMatrix >

    where 
        StrictlyUpperTriangularMatrix:      MatrixAlgebra< 
                                                RowIndex        =   Index,          // we use the Index type parameter to ensure that row and column indices to have the same type
                                                ColumnIndex     =   Index,          // we use the Index type parameter to ensure that row and column indices to have the same type
                                                RowEntry:           KeyValNew,
                                                ColumnEntry:        KeyValNew
                                            >, 
        Index:                              Clone + Eq + std::fmt::Debug,

{
    type RingOperator                       =   StrictlyUpperTriangularMatrix::RingOperator                  ;

    type OrderOperatorForRowEntries         =   StrictlyUpperTriangularMatrix::OrderOperatorForRowEntries    ;

    type OrderOperatorForRowIndices         =   StrictlyUpperTriangularMatrix::OrderOperatorForRowIndices    ;

    type OrderOperatorForColumnEntries      =   StrictlyUpperTriangularMatrix::OrderOperatorForColumnEntries ;

    type OrderOperatorForColumnIndices      =   StrictlyUpperTriangularMatrix::OrderOperatorForColumnIndices ;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.strictly_upper_triangular_matrix.ring_operator()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.strictly_upper_triangular_matrix.order_operator_for_row_entries()
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.strictly_upper_triangular_matrix.order_operator_for_row_indices()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.strictly_upper_triangular_matrix.order_operator_for_column_entries()
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.strictly_upper_triangular_matrix.order_operator_for_column_indices()
    }
}






//  Matrix Oracle Operations
//  ---------------------------------------------------------------------------


impl < StrictlyUpperTriangularMatrix > 

    MatrixOracleOperations for 

    SumOfScalarAndStrictlyUpperTriangularMatrices 
        < StrictlyUpperTriangularMatrix >

    where 
        StrictlyUpperTriangularMatrix:          MatrixOracle // required in type definition        
{}














//  TESTS
//  =========================================================================================================

//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    use std::iter::Cloned;

    

    


    
    #[test] 
    fn test_PrependEntryToRow() {

        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::types::scalar_diagonal_triangle::SumOfScalarAndStrictlyUpperTriangularMatrices;
        use crate::algebra::matrices::query::MatrixOracle;
        use std::slice::Iter;

        type Row<'a> = Cloned< Iter< 'a, (i32, i32) > >;
        let strictly_upper_triangular_matrix  =   VecOfVec::new(
                                            vec![
                                                vec![ (1, 1), (2, 1) ],
                                                vec![         (2, 1) ],
                                                vec![                ],
                                            ]
                                        ).ok().unwrap();

        let matrix_prepended
        // :   SumOfScalarAndStrictlyUpperTriangularMatrices< 
        //                                 & VecOfVec<i32, i32>, 
        //                                 Cloned< Iter< (i32, i32) > >, 
        //                                 i32, 
        //                                 i32
        //                             >
            =   SumOfScalarAndStrictlyUpperTriangularMatrices::new( & strictly_upper_triangular_matrix, 1 );

        let intended_result  =   VecOfVec::new(
                                            vec![
                                                vec![ (0, 1), (1, 1), (2, 1) ],
                                                vec![         (1, 1), (2, 1) ],
                                                vec![                 (2, 1) ],
                                            ]
                                        ).ok().unwrap();     
        // let intended_result_ref: &'a VecOfVec< i32, i32 > = & &intended_result;
                                        
        for row_index in 0 .. 3 {
            itertools::assert_equal( 
                matrix_prepended.row( &row_index ),
                (& intended_result).row( &row_index ),
            )
        }                                              

    }

}