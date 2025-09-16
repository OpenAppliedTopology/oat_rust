//! Wrapper for several pieces of data commonly used together with a matrix

use derive_getters::Dissolve;
use num::rational::Ratio;

use crate::{algebra::{matrices::{operations::MatrixOracleOperations, query::{ MatrixAlgebra, MatrixOracle,}}, rings::{traits::SemiringOperations, types::{field_prime_order::BooleanField, native::{FieldFloat64, FieldRationalSize, RingOperatorForNativeRustNumberType}}}}, utilities::order::{JudgeOrder, OrderOperatorAuto, OrderOperatorByKey}};





/// Wrapper for several pieces of data commonly used together with a matrix
/// 
/// The user-provided matrix must be compatible with the user-provided ring and order operators, as
/// per the guidelines provided for the [MatrixAlgebra] trait. 
/// **If these conditions are violated then any calculations performed with the [MatrixAlgebraPacket] may be incorrect.**
/// **Moreover, it's possible that no errors or warnings will be generated. Therefore use this object with caution.**
#[derive(Clone,Copy,Debug,Dissolve,Eq,PartialEq,Ord,PartialOrd)]
pub struct MatrixAlgebraPacket
                < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, > 
{
    pub matrix:                             Matrix,
    pub ring_operator:                      RingOperator,
    pub order_operator_for_row_entries:     OrderOperatorForRowEntries,
    pub order_operator_for_row_indices:     OrderOperatorForRowIndices,    
    pub order_operator_for_column_entries:  OrderOperatorForColumnEntries,
    pub order_operator_for_column_indices:  OrderOperatorForColumnIndices,        
}

impl < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, >

    MatrixAlgebraPacket
        < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, >

     
{
    /// A reference to the matrix 
    pub fn matrix_ref(&self) -> & Matrix { & self.matrix }
}




//  SPECIALIZED IMPLEMENTATION FOR SIMPLE TYPES


impl < Matrix, RingOperator, >

    MatrixAlgebraPacket
        < 
            Matrix, 
            RingOperator, 
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
        >
{
    /// Wraps `matrix` in a [`MatrixAlgebraPacket`]
    /// 
    /// This method is only available for matrices whose row and column indices implement [PartialOrd]. The constructor will assign the order from
    /// `PartialOrd` to the row and column entries and indices. The user specifies the ring operator.
    pub fn with_default_order( matrix: Matrix, ring_operator: RingOperator ) -> Self { 
        MatrixAlgebraPacket{
            matrix,
            ring_operator,
            order_operator_for_row_entries: OrderOperatorByKey::new(),
            order_operator_for_row_indices: OrderOperatorAuto,
            order_operator_for_column_entries: OrderOperatorByKey::new(),
            order_operator_for_column_indices: OrderOperatorAuto,          
        }
    }
}


impl < Matrix >

    MatrixAlgebraPacket
        < 
            Matrix, 
            BooleanField, 
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
        >
    where
        Matrix:         MatrixOracle< 
                            Coefficient         =   bool,
                            RowIndex:               PartialOrd,
                            ColumnIndex:            PartialOrd,
                        >,   
{
    /// Wraps the input `matrix`  with boolean coefficients in a [`MatrixAlgebraPacket`]
    /// 
    /// This method is only available for matrices whose row and column indices implement [PartialOrd]. The constructor will assign the order from
    /// `PartialOrd` to the row and column entries and indices. The packets' ring operator will be [`BooleanField`]    
    pub fn with_default_order_and_boolean_coefficients( matrix: Matrix ) -> Self { 
        MatrixAlgebraPacket{
            matrix,
            ring_operator: BooleanField::new(),
            order_operator_for_row_entries: OrderOperatorByKey::new(),
            order_operator_for_row_indices: OrderOperatorAuto,
            order_operator_for_column_entries: OrderOperatorByKey::new(),
            order_operator_for_column_indices: OrderOperatorAuto,          
        }
    }
}




impl < Matrix >

    MatrixAlgebraPacket
        < 
            Matrix, 
            RingOperatorForNativeRustNumberType<i64>, 
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
        >
    where
        Matrix:         MatrixOracle< 
                            Coefficient         =   i64,
                            RowIndex:               PartialOrd,
                            ColumnIndex:            PartialOrd,
                        >,   
{
    /// Wraps the input `matrix`  with integer `i64` coefficients in a [`MatrixAlgebraPacket`]
    /// 
    /// This method is only available for matrices whose row and column indices implement [PartialOrd]. The constructor will assign the order from
    /// `PartialOrd` to the row and column entries and indices. The packets' ring operator will be [`BooleanField`]    
    pub fn with_default_order_and_i64_coefficients( matrix: Matrix ) -> Self { 
        MatrixAlgebraPacket{
            matrix,
            ring_operator: RingOperatorForNativeRustNumberType::<i64>::new(),
            order_operator_for_row_entries: OrderOperatorByKey::new(),
            order_operator_for_row_indices: OrderOperatorAuto,
            order_operator_for_column_entries: OrderOperatorByKey::new(),
            order_operator_for_column_indices: OrderOperatorAuto,          
        }
    }
}











impl < Matrix >

    MatrixAlgebraPacket
        < 
            Matrix, 
            RingOperatorForNativeRustNumberType<usize>, 
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
        >
    where
        Matrix:         MatrixOracle< 
                            Coefficient         =   i64,
                            RowIndex:               PartialOrd,
                            ColumnIndex:            PartialOrd,
                        >,   
{
    /// Wraps the input `matrix`  with integer `usize` coefficients in a [`MatrixAlgebraPacket`]
    /// 
    /// This method is only available for matrices whose row and column indices implement [PartialOrd]. The constructor will assign the order from
    /// `PartialOrd` to the row and column entries and indices. The packets' ring operator will be [`BooleanField`]    
    pub fn with_default_order_and_usize_coefficients( matrix: Matrix ) -> Self { 
        MatrixAlgebraPacket{
            matrix,
            ring_operator: RingOperatorForNativeRustNumberType::<usize>::new(),
            order_operator_for_row_entries: OrderOperatorByKey::new(),
            order_operator_for_row_indices: OrderOperatorAuto,
            order_operator_for_column_entries: OrderOperatorByKey::new(),
            order_operator_for_column_indices: OrderOperatorAuto,          
        }
    }
}













impl < Matrix, >

    MatrixAlgebraPacket
        < 
            Matrix, 
            FieldRationalSize, 
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
        >
    where
        Matrix:         MatrixOracle< 
                            Coefficient         =   Ratio< isize >,
                            RowIndex:               PartialOrd,
                            ColumnIndex:            PartialOrd,
                        >,                  
{
    /// Wraps `matrix`  with rational coefficients in a [`MatrixAlgebraPacket`]
    /// 
    /// This method is only available for matrices whose row and column indices implement [PartialOrd]. The constructor will assign the order from
    /// `PartialOrd` to the row and column entries and indices. The packets' ring operator will be [`FieldRationalSize`]
    pub fn with_default_order_and_rational_coefficients( matrix: Matrix ) -> Self { 
        MatrixAlgebraPacket{
            matrix,
            ring_operator: FieldRationalSize::new(),
            order_operator_for_row_entries: OrderOperatorByKey::new(),
            order_operator_for_row_indices: OrderOperatorAuto,
            order_operator_for_column_entries: OrderOperatorByKey::new(),
            order_operator_for_column_indices: OrderOperatorAuto,          
        }
    }
}


impl < Matrix, >

    MatrixAlgebraPacket
        < 
            Matrix, 
            FieldFloat64, 
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
            OrderOperatorByKey, // automatically compares entries by key
            OrderOperatorAuto,  // compares indices using the default order
        >
    where
        Matrix:         MatrixOracle< 
                            Coefficient         =   f64,
                            RowIndex:               PartialOrd,
                            ColumnIndex:            PartialOrd,
                        >,                     
{
    /// Wraps `matrix`  with rational coefficients in a [`MatrixAlgebraPacket`]
    /// 
    /// This method is only available for matrices whose row and column indices implement [PartialOrd]. The constructor will assign the order from
    /// `PartialOrd` to the row and column entries and indices. The packets' ring operator will be [`FieldRationalSize`]
    pub fn with_default_order_and_f64_coefficients( matrix: Matrix ) -> Self { 
        MatrixAlgebraPacket{
            matrix,
            ring_operator: FieldFloat64::new(),
            order_operator_for_row_entries: OrderOperatorByKey::new(),
            order_operator_for_row_indices: OrderOperatorAuto,
            order_operator_for_column_entries: OrderOperatorByKey::new(),
            order_operator_for_column_indices: OrderOperatorAuto,          
        }
    }
}










impl < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, > 

    MatrixOracle for

    MatrixAlgebraPacket
        < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, > 

    where 
        Matrix:                         MatrixOracle,

{
    
    type Coefficient            =   Matrix::Coefficient  ;    // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   Matrix::RowIndex     ;    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   Matrix::ColumnIndex  ;    // The type of column indices    
    type RowEntry               =   Matrix::RowEntry     ;  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   Matrix::ColumnEntry  ;  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =   Matrix::Row          ;  // What you get when you ask for a row.
    type RowReverse             =   Matrix::RowReverse   ;  // What you get when you ask for a row with the order of entries reversed
    type Column                 =   Matrix::Column       ;  // What you get when you ask for a column
    type ColumnReverse          =   Matrix::ColumnReverse;  // What you get when you ask for a column with the order of entries reversed                             

    /// Should return `Some(x)` as long as the row and column indices are valid.
    fn structural_nonzero_entry(                   &   self, row:   &Self::RowIndex, column: &Self::ColumnIndex ) ->  Option< Self::Coefficient >  { self.matrix.structural_nonzero_entry(row,column) }
    fn has_column_for_index(  &   self, index: &Self::ColumnIndex)    -> bool                         { self.matrix.has_column_for_index(index) }
    fn has_row_for_index(     &   self, index: &Self::RowIndex)       -> bool                         { self.matrix.has_row_for_index(index) }    
    fn row(                     &   self, index: &Self::RowIndex   )    -> Self::Row                    { self.matrix.row( index )                }
    fn row_result(                 &   self, index: &Self::RowIndex   )    -> Result< Self::Row, Self::RowIndex >            { self.matrix.row_result( index )            }
    fn row_reverse(             &   self, index: &Self::RowIndex   )    -> Self::RowReverse             { self.matrix.row_reverse( index )        }
    fn row_reverse_result(         &   self, index: &Self::RowIndex   )    -> Result< Self::RowReverse, Self::RowIndex >     { self.matrix.row_reverse_result( index )    }
    fn column(                  &   self, index: &Self::ColumnIndex)    -> Self::Column                 { self.matrix.column( index )             }
    fn column_result(              &   self, index: &Self::ColumnIndex)    -> Result< Self::Column, Self::ColumnIndex >         { self.matrix.column_result( index )         }
    fn column_reverse(          &   self, index: &Self::ColumnIndex)    -> Self::ColumnReverse          { self.matrix.column_reverse( index )     }
    fn column_reverse_result(      &   self, index: &Self::ColumnIndex)    -> Result< Self::ColumnReverse, Self::ColumnIndex >  { self.matrix.column_reverse_result( index ) }

} 





impl < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, > 

    MatrixAlgebra for

    MatrixAlgebraPacket
        < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, > 

    where 
        Matrix:                         MatrixOracle,
        RingOperator:                   Clone + SemiringOperations< Element = Self::Coefficient >,
        OrderOperatorForRowEntries:     Clone + JudgeOrder< Self::RowEntry >,
        OrderOperatorForRowIndices:     Clone + JudgeOrder< Self::RowIndex >,
        OrderOperatorForColumnEntries:  Clone + JudgeOrder< Self::ColumnEntry >,
        OrderOperatorForColumnIndices:  Clone + JudgeOrder< Self::ColumnIndex >,   
    
{
    type RingOperator                   =   RingOperator                    ;
    type OrderOperatorForRowEntries     =   OrderOperatorForRowEntries      ;
    type OrderOperatorForRowIndices     =   OrderOperatorForRowIndices      ;
    type OrderOperatorForColumnEntries  =   OrderOperatorForColumnEntries   ;
    type OrderOperatorForColumnIndices  =   OrderOperatorForColumnIndices   ;

    fn ring_operator( &self ) -> Self::RingOperator {                                           self.ring_operator.clone()                      }
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {            self.order_operator_for_row_entries.clone()     }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {            self.order_operator_for_row_indices.clone()     }    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {      self.order_operator_for_column_entries.clone()  }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {      self.order_operator_for_column_indices.clone()  }    
}







impl < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, > 

    MatrixOracleOperations for

    MatrixAlgebraPacket
        < Matrix, RingOperator, OrderOperatorForRowEntries, OrderOperatorForRowIndices, OrderOperatorForColumnEntries, OrderOperatorForColumnIndices, > 
{}        