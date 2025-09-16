use std::iter::Peekable;

use derive_new::new;
use derive_getters::{Dissolve, Getters};

use crate::algebra::matrices::{operations::MatrixOracleOperations, query::{MatrixAlgebra, MatrixOracle}};





/// Wraps the rows and columns of a matrix in a [Peekable] struct, so that they implement [PeekUnqualified](crate::utilities::iterators::general::PeekUnqualified)
/// 
/// A [Peekable] is simply a wrapper for an iterator `I` that holds both `I` and one of the items of `I`.
/// This allows the user to "peek" at the next item which `I` should return, without removing that element.
/// NB: the [Peekable] object has a native method called `peek`, but there is a trait called [PeekUnqualified](crate::utilities::iterators::general::PeekUnqualified)
/// in OAT which provides an equivalent method called `peek_unqualified`; unlike `peek`, which is type-specific,
/// `peek_unqualified` can be implemented on a variety of types.
/// 
/// # Why does this struct exist?
/// 
/// This struct incurs a minor memory expense, because it stores an item in addition to the iterator. However,
/// many useful operations require the ability to peek at the next element of an iterator. 
#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd,new)]
pub struct PeekableMatrix< T > { 
    pre_peekable_matrix: T 
}



impl < T >

    MatrixOracle for

    PeekableMatrix < T >

    where
        T:          MatrixOracle

{
    type Coefficient            =   T::Coefficient  ;
    type RowIndex               =   T::RowIndex     ;
    type ColumnIndex            =   T::ColumnIndex  ;
    type RowEntry               =   T::RowEntry     ;
    type ColumnEntry            =   T::ColumnEntry  ;
    type Row                    =   Peekable< T::Row           >;
    type RowReverse             =   Peekable< T::RowReverse    >;
    type Column                 =   Peekable< T::Column        >;
    type ColumnReverse          =   Peekable< T::ColumnReverse >;



    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        self.pre_peekable_matrix.row( index ).peekable()
    }

    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        self.pre_peekable_matrix.row_reverse( index ).peekable()
    }

    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
        self.pre_peekable_matrix.column( index ).peekable()
    }

    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        self.pre_peekable_matrix.column_reverse( index ).peekable()
    }

    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.pre_peekable_matrix.has_row_for_index(index)
    }

    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.pre_peekable_matrix.has_column_for_index(index)
    }

    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        self.pre_peekable_matrix.structural_nonzero_entry( row, column )
    }
}    




impl < T >

    MatrixAlgebra for

    PeekableMatrix < T >

    where
        T:          MatrixAlgebra

{
    type RingOperator                       =   T::RingOperator                 ;
    type OrderOperatorForRowEntries         =   T::OrderOperatorForRowEntries   ;
    type OrderOperatorForRowIndices         =   T::OrderOperatorForRowIndices   ;
    type OrderOperatorForColumnEntries      =   T::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices      =   T::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.pre_peekable_matrix.ring_operator()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        self.pre_peekable_matrix.order_operator_for_row_entries()
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        self.pre_peekable_matrix.order_operator_for_row_indices()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        self.pre_peekable_matrix.order_operator_for_column_entries()
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        self.pre_peekable_matrix.order_operator_for_column_indices()
    }
}    









impl < T >

    MatrixOracleOperations for

    PeekableMatrix < T >
{}    