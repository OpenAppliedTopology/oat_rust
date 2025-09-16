//! Lazy transpose and anti-transpose; wraps around another matrix, and swaps order and/or major vs. columns.

use derive_getters::{Getters, Dissolve};
use derive_new::new;

use crate::algebra::matrices::query::{MatrixOracle, MatrixAlgebra};
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::chain_complexes::ChainComplex;
use crate::utilities::order::ReverseOrder;






//  ANTITRANSPOSE
//  -----------------------------------------------------------------------------------------------

/// Order anti-tranpose of a matrix, evaluated in a lazy fashion. (caution: this doesn't work quite the same as normal anti-transpose for matrices indexed by integers)
/// 
/// Concretely, if `m` is a matrix and `a = OrderAntiTranspose::new(m)` then
/// - `a.row(&index)`         returns `m.column_reverse(&index)`
/// - `a.row_reverse(&index)` returns `m.column(&index)`
/// - `a.column(&index)` returns         `m.row_reverse(&index)`
/// - `a.column_reverse(&index)` returns `m.row(&index)`
/// - `a.structural_nonzero_entry(&i, &j)` returns `m.structural_nonzero_entry(&j, &i)`
/// 
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  
/// 
/// 
/// # Caution
/// 
/// This operation doesn't work quite the same as for normal anti-transpose for matrices indexed by integers.
/// **The key difference is that [OrderAntiTranspose] does not change the index of any nonzero entry. All it does is change the order in which row and column iterators return entries.**
/// 
/// To illustrate, let `m` be the following matrix:
/// 
/// ```text
/// 0 5
/// 0 0
/// 0 6
/// 7 8
/// ```
/// 
/// Let `x` be the normal anti-transpose of `m`, which is 
/// 
/// ```text
/// 8 6 0 5
/// 7 0 0 0
/// ```
/// 
/// Finally, suppose that `o = OrderAntiTranspose::new(m)`. Then
/// 
/// - `o.structural_nonzero_entry(&0,&0) = None` because `o.structural_nonzero_entry(&0,&0) = m.structural_nonzero_entry(&0,&0)`. However,
/// - `x.structural_nonzero_entry(&0,&0) = Some(8)`
/// - `o.row(&1) = [(3,8),(2,6),(0,5)]`, however
/// - `x.row(&1) = [(0,7)]`
/// 
/// **If you need to obtain `x` instead of `o`, see the [antitranspose_deep](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep) method
/// for [VecOfVec](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec) datastructures.** Note that this operation is only available 
/// for matrices with rows and columns indexed by integers.
/// 
/// 
/// 
/// # Matrix operators
/// 
/// If `Matrix` implements [MatrixAlgebra] then so does `OrderAntiTranspose<Matrix>`. 
/// - `OrderAntiTranspose<Matrix>` has the same ring operator as `Matrix`
/// - `OrderAntiTranspose<Matrix>` has all the same order operators as `Matrix`, but *reversed*. For example, if `Matrix::OrderOperatorForRowEntries` has type `T`, then
/// `OrderAntiTranspose<Matrix>::OrderOperatorForRowOperators` has type `ReverseOrder<T>`.
/// 
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::transpose::{Transpose, OrderAntiTranspose};
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::utilities::iterators::general::result_iterators_are_elementwise_equal;
/// 
/// // matrix
/// let matrix: & VecOfVec<usize,usize> =   & VecOfVec::new( vec![
///                                             vec![ (0,0), (1,1), (2,2) ],
///                                             vec![ (0,3), (1,4), (2,5) ],
///                                         ] ).ok().unwrap();
/// 
/// // its transpose
/// let tran    =   Transpose::new( &matrix );
/// 
/// // its antitranspose
/// let anti    =   OrderAntiTranspose::new( &matrix );
/// 
/// 
/// // check that rows and columns are (anti)transposed correctly
/// for index in 0 .. 2 {
///     assert!( matrix.row(&index).eq(                  tran.column(&index) )               );
///     assert!( matrix.row_reverse(&index).eq(          tran.column_reverse(&index) )       );            
///     assert!( matrix.row(&index).eq(                  anti.column_reverse(&index) )       );
///     assert!( matrix.row_reverse(&index).eq(          anti.column(&index) )               );
/// }
/// 
/// for index in 0 .. 3 {
///     assert!( matrix.column(&index).eq(                tran.row(&index) )                 );
///     assert!( matrix.column_reverse(&index).eq(        tran.row_reverse(&index) )         );
///     assert!( matrix.column(&index).eq(                anti.row_reverse(&index) )         );            
///     assert!( matrix.column_reverse(&index).eq(        anti.row(&index) )                 );
/// }        
///                                                                                             
/// // check that rows and columns are (anti)transposed correctly, even when we pass invalid indices     
/// for index in 0 .. 5 {    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_result(&index),                 tran.column_result(&index) 
///     ));    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_result(&index),                 anti.column_reverse_result(&index) 
///     ));    
/// 
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_reverse_result(&index),         tran.column_reverse_result(&index) 
///     ));  
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_reverse_result(&index),         anti.column_result(&index) 
///     ));                                       
/// }
/// 
/// for index in 0 .. 5 {    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_result(&index),                 tran.row_result(&index) 
///     ));    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_result(&index),                 anti.row_reverse_result(&index) 
///     ));    
/// 
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_reverse_result(&index),         tran.row_reverse_result(&index) 
///     ));  
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_reverse_result(&index),         anti.row_result(&index) 
///     ));                                       
/// }
/// ```
#[derive(Debug, Copy, Clone, Eq, PartialEq, Getters, Dissolve, new)]
pub struct OrderAntiTranspose< Matrix > { matrix_to_antitranspose: Matrix }



//  MATRIX ORACLE
//  ---------------------------

impl< Matrix > 

    MatrixOracle for 
    
    OrderAntiTranspose< Matrix >

    where
        Matrix:     MatrixOracle
{
    type Coefficient        =   Matrix::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   Matrix::ColumnIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   Matrix::RowIndex;       // The type of column indices
    
    type RowEntry           =   Matrix::ColumnEntry;          // The type of entries in each row; these are essentially pairs of form `(&column_index, coefficient)`
    type ColumnEntry        =   Matrix::RowEntry;       // The type of entries in each column; these are essentially pairs of form `(&row_index, coefficient)`
    
    type Row                =   Matrix::ColumnReverse;               // What you get when you ask for a row.
    type RowReverse         =   Matrix::Column;        // What you get when you ask for a row with the order of entries reversed
    
    type Column             =   Matrix::RowReverse;            // What you get when you ask for a column
    type ColumnReverse      =   Matrix::Row;     // What you get when you ask for a column with the order of entries reversed 

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        self.matrix_to_antitranspose.structural_nonzero_entry( column, row )
    }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool 
        { self.matrix_to_antitranspose.has_column_for_index(index) }        
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool 
        { self.matrix_to_antitranspose.has_row_for_index(index) }    

    fn row(                     &self,  index: &Self::RowIndex    )       -> Self::Row 
        { self.matrix_to_antitranspose.column_reverse(index) }
    fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex >
        { self.matrix_to_antitranspose.column_reverse_result(index) }     
    fn row_reverse(             &self,  index: &Self::RowIndex    )       -> Self::RowReverse
        { self.matrix_to_antitranspose.column(index) }    
    fn row_reverse_result(         &   self, index: & Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex >
        { self.matrix_to_antitranspose.column_result(index) }
    
    fn column(                  &self,  index: &Self::ColumnIndex )       -> Self::Column
        { self.matrix_to_antitranspose.row_reverse(index) }
    fn column_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex >
        { self.matrix_to_antitranspose.row_reverse_result(index) }    
    fn column_reverse(          &self,  index: &Self::ColumnIndex )       -> Self::ColumnReverse
        { self.matrix_to_antitranspose.row(index) }
    fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex >
        { self.matrix_to_antitranspose.row_result(index) }
} 


//  MATRIX ALGEBRA
//  --------------------------


impl < Matrix >

    MatrixAlgebra for

    OrderAntiTranspose
        < Matrix > where

    Matrix: MatrixAlgebra 
{
    type OrderOperatorForColumnEntries =   ReverseOrder< Matrix::OrderOperatorForRowEntries >;
    type OrderOperatorForColumnIndices =   ReverseOrder< Matrix::OrderOperatorForRowIndices >;
    type OrderOperatorForRowEntries    =   ReverseOrder< Matrix::OrderOperatorForColumnEntries >;
    type OrderOperatorForRowIndices    =   ReverseOrder< Matrix::OrderOperatorForColumnIndices >;
    type RingOperator               =   Matrix::RingOperator;

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { 
        ReverseOrder::new( self.matrix_to_antitranspose.order_operator_for_row_entries() )
    }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { 
        ReverseOrder::new( self.matrix_to_antitranspose.order_operator_for_row_indices() )
    }   
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { 
        ReverseOrder::new( self.matrix_to_antitranspose.order_operator_for_column_entries() )
    }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { 
        ReverseOrder::new( self.matrix_to_antitranspose.order_operator_for_column_indices() )
    }
    fn ring_operator( &self ) -> Self::RingOperator {
        self.matrix_to_antitranspose.ring_operator()
    }     
}



//  MATRIX ORACLE OPERATIONS
//  ---------------------------

impl< Matrix > 

    MatrixOracleOperations for 
    
    OrderAntiTranspose< Matrix >
{}






//  CHAIN COMPLEX
//  --------------------------


impl < Matrix, IndexForRowsAndColumns >

    ChainComplex for

    OrderAntiTranspose
        < Matrix > where

    Matrix:     ChainComplex + MatrixOracle< 
                    RowIndex =  IndexForRowsAndColumns, 
                    ColumnIndex = IndexForRowsAndColumns,
                >,
{

    type BasisVectorIndicesIterable = Vec< Self::RowIndex >; 

    /// Returns an iterable that runs over the indices of the basis vectors.
    /// 
    /// Indices should be returned sorted in the same order they are asigned in the rows/columns of the differential matrix.
    /// That is, strictly increasing order with respect to `self.order_operator_for_row_indices()`.
    fn basis_vector_indices_for_dimension( &self, dimension: isize ) -> Self::BasisVectorIndicesIterable {
        let mut indices: Vec<_> = self.matrix_to_antitranspose.basis_vector_indices_for_dimension( dimension )
            .into_iter()
            .collect();
        indices.reverse();
        indices
    }


    /// Returns the dimension of the basis vector with the given index.
    fn dimension_for_basis_vector_with_index( &self, index: & Self::RowIndex ) -> Result<isize, Self::RowIndex> {
        self.matrix_to_antitranspose.dimension_for_basis_vector_with_index( index )
    }


}







//  TRANSPOSE
//  -----------------------------------------------------------------------------------------------

/// Transpose of a matrix, evaluated in a lazy fashion.
/// 
/// This struct is "lazy" in the sense that it does not transform the underlying data; rather, it simply swaps access commands with other
/// access commands.  For example, to execute the method `row`, a `Transpose` struct simply calls `column` on
/// the underlying matrix.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::matrices::types::transpose::{Transpose, OrderAntiTranspose};
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::query::MatrixOracle;
/// use oat_rust::utilities::iterators::general::result_iterators_are_elementwise_equal;
/// 
/// // matrix
/// let matrix: & VecOfVec<usize,usize> =   & VecOfVec::new( vec![
///                                             vec![ (0,0), (1,1), (2,2) ],
///                                             vec![ (0,3), (1,4), (2,5) ],
///                                         ] ).ok().unwrap();
/// 
/// // its transpose
/// let tran    =   Transpose::new( &matrix );
/// 
/// // its antitranspose
/// let anti    =   OrderAntiTranspose::new( &matrix );
/// 
/// 
/// // check that rows and columns are (anti)transposed correctly
/// for index in 0 .. 2 {
///     assert!( matrix.row(&index).eq(                  tran.column(&index) )               );
///     assert!( matrix.row_reverse(&index).eq(          tran.column_reverse(&index) )       );            
///     assert!( matrix.row(&index).eq(                  anti.column_reverse(&index) )       );
///     assert!( matrix.row_reverse(&index).eq(          anti.column(&index) )               );
/// }
/// 
/// for index in 0 .. 3 {
///     assert!( matrix.column(&index).eq(                tran.row(&index) )                 );
///     assert!( matrix.column_reverse(&index).eq(        tran.row_reverse(&index) )         );
///     assert!( matrix.column(&index).eq(                anti.row_reverse(&index) )         );            
///     assert!( matrix.column_reverse(&index).eq(        anti.row(&index) )                 );
/// }        
///                                                                                             
/// // check that rows and columns are (anti)transposed correctly, even when we pass invalid indices     
/// for index in 0 .. 5 {    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_result(&index),                 tran.column_result(&index) 
///     ));    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_result(&index),                 anti.column_reverse_result(&index) 
///     ));    
/// 
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_reverse_result(&index),         tran.column_reverse_result(&index) 
///     ));  
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.row_reverse_result(&index),         anti.column_result(&index) 
///     ));                                       
/// }
/// 
/// for index in 0 .. 5 {    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_result(&index),                 tran.row_result(&index) 
///     ));    
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_result(&index),                 anti.row_reverse_result(&index) 
///     ));    
/// 
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_reverse_result(&index),         tran.row_reverse_result(&index) 
///     ));  
///     assert!( result_iterators_are_elementwise_equal( 
///         matrix.column_reverse_result(&index),         anti.row_result(&index) 
///     ));                                       
/// }
/// ```
/// 
#[derive(Debug, Copy, Clone, Eq, PartialEq, Getters, Dissolve, new)]
pub struct Transpose< Matrix > { untransposed: Matrix }


//  MATRIX ORACLE
//  --------------------------


impl< Matrix > 

MatrixOracle for 

Transpose< Matrix >

where
    Matrix:     MatrixOracle

{
    type Coefficient        =   Matrix::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex           =   Matrix::ColumnIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex        =   Matrix::RowIndex;       // The type of column indices
    
    type RowEntry           =   Matrix::ColumnEntry;          // The type of entries in each row; these are essentially pairs of form `(&column_index, coefficient)`
    type ColumnEntry        =   Matrix::RowEntry;       // The type of entries in each column; these are essentially pairs of form `(&row_index, coefficient)`
    
    type Row                =   Matrix::Column;               // What you get when you ask for a row.
    type RowReverse         =   Matrix::ColumnReverse;        // What you get when you ask for a row with the order of entries reversed
    
    type Column             =   Matrix::Row;            // What you get when you ask for a column
    type ColumnReverse      =   Matrix::RowReverse;     // What you get when you ask for a column with the order of entries reversed 

    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > 
        { self.untransposed.structural_nonzero_entry( column, row ) }
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool 
        { self.untransposed.has_row_for_index(index) }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool 
        { self.untransposed.has_column_for_index(index) }

    fn row(                     &self,  index: &Self::RowIndex    )       -> Self::Row 
        { self.untransposed.column(index) }
    fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex >
        { self.untransposed.column_result(index) }     
    fn row_reverse(             &self,  index: &Self::RowIndex    )       -> Self::RowReverse
        { self.untransposed.column_reverse(index) }    
    fn row_reverse_result(         &   self, index: & Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex >
        { self.untransposed.column_reverse_result(index) }
    
    fn column(                  &self,  index: &Self::ColumnIndex )       -> Self::Column
        { self.untransposed.row(index) }
    fn column_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex >
        { self.untransposed.row_result(index) }    
    fn column_reverse(          &self,  index: &Self::ColumnIndex )       -> Self::ColumnReverse
        { self.untransposed.row_reverse(index) }
    fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex >
        { self.untransposed.row_reverse_result(index) }
} 



//  MATRIX ALGEBRA
//  --------------------------



impl < Matrix >

    MatrixAlgebra for

    Transpose< Matrix > where

    Matrix: MatrixAlgebra 
{
    type OrderOperatorForColumnEntries =   Matrix::OrderOperatorForRowEntries;
    type OrderOperatorForColumnIndices =   Matrix::OrderOperatorForRowIndices;
    type OrderOperatorForRowEntries    =   Matrix::OrderOperatorForColumnEntries;
    type OrderOperatorForRowIndices    =   Matrix::OrderOperatorForColumnIndices;
    type RingOperator               =   Matrix::RingOperator;

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { 
        self.untransposed.order_operator_for_row_entries()
    }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { 
        self.untransposed.order_operator_for_row_indices()
    }   
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { 
        self.untransposed.order_operator_for_column_entries()
    }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { 
        self.untransposed.order_operator_for_column_indices()
    }
    fn ring_operator( &self ) -> Self::RingOperator {
        self.untransposed.ring_operator()
    }     
}





//  MATRIX ORACLE OPERATIONS
//  --------------------------


impl< Matrix > 

MatrixOracleOperations for 

Transpose< Matrix >
{}








//  ---------------------------------------------------------------------------
//  TEST
//  ---------------------------------------------------------------------------


#[cfg(test)]
mod tests {
    

    


    #[test]
    fn test_transpose_and_antitranspose() {
        
        use crate::algebra::matrices::types::transpose::{Transpose, OrderAntiTranspose};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        use itertools;
        
        // matrix
        let a     =   VecOfVec::new( vec![
                                                vec![ (0,0), (1,1), (2,2) ],
                                                vec![ (0,3), (1,4), (2,5) ],
                                            ] ).ok().unwrap();
        
        // its transpose
        let tran    =   Transpose::new( &a );

        // its antitranspose
        let anti    =   OrderAntiTranspose::new( &a );                                                                                                                                 
        

        for row in 0 .. 2 {
            
            assert!( itertools::equal( (&a).row(&row),                           tran.column(&row ) )                        );
            assert!( itertools::equal( (&a).row_result(&row).unwrap(),              tran.column_result(&row ).unwrap() )           );
            assert!( itertools::equal( (&a).row_reverse(&row),                   tran.column_reverse(&row))                  );
            assert!( itertools::equal( (&a).row_reverse_result(&row).unwrap(),      tran.column_reverse_result(&row).unwrap() )    );

            assert!( itertools::equal( (&a).row_reverse(&row),                   anti.column(&row) )                         );
            assert!( itertools::equal( (&a).row_reverse_result(&row).unwrap(),      anti.column_result(&row).unwrap() )            ); 
            assert!( itertools::equal( (&a).row(&row),                           anti.column_reverse(&row) )                 );
            assert!( itertools::equal( (&a).row_result(&row).unwrap(),              anti.column_reverse_result(&row).unwrap() )    );
        }

        for col in 0 .. 3 {

            assert!( itertools::equal( (&a).column(&col),                        tran.row(&col) )                             );
            assert!( itertools::equal( (&a).column_result(&col).unwrap(),           tran.row_result(&col).unwrap() )                );
            assert!( itertools::equal( (&a).column_reverse(&col),                tran.row_reverse(&col) )                     );
            assert!( itertools::equal( (&a).column_reverse_result(&col).unwrap(),   tran.row_reverse_result(&col).unwrap() )        );
            
            assert!( itertools::equal( (&a).column(&col),                        anti.row_reverse(&col) )                     );
            assert!( itertools::equal( (&a).column_result(&col).unwrap(),           anti.row_reverse_result(&col).unwrap() )        );
            assert!( itertools::equal( (&a).column_reverse(&col),                anti.row(&col) )                             );
            assert!( itertools::equal( (&a).column_reverse_result(&col).unwrap(),   anti.row_result(&col).unwrap() )                );
            
        }    

    }    


    #[test]
    fn test_antitranspose() {
        
        use crate::algebra::matrices::types::transpose::{Transpose, OrderAntiTranspose};
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::utilities::iterators::general::result_iterators_are_elementwise_equal;

        // matrix
        let matrix: & VecOfVec<usize,usize> =   & VecOfVec::new( vec![
                                                    vec![ (0,0), (1,1), (2,2) ],
                                                    vec![ (0,3), (1,4), (2,5) ],
                                                ] ).ok().unwrap();
        
        // its transpose
        let tran    =   Transpose::new( &matrix );
        
        // its antitranspose
        let anti    =   OrderAntiTranspose::new( &matrix );


        // check that rows and columns are (anti)transposed correctly
        for index in 0 .. 2 {
            assert!( matrix.row(&index).eq(                  tran.column(&index) )               );
            assert!( matrix.row_reverse(&index).eq(          tran.column_reverse(&index) )       );            
            assert!( matrix.row(&index).eq(                  anti.column_reverse(&index) )       );
            assert!( matrix.row_reverse(&index).eq(          anti.column(&index) )               );
        }

        for index in 0 .. 3 {
            assert!( matrix.column(&index).eq(                tran.row(&index) )                 );
            assert!( matrix.column_reverse(&index).eq(        tran.row_reverse(&index) )         );
            assert!( matrix.column(&index).eq(                anti.row_reverse(&index) )         );            
            assert!( matrix.column_reverse(&index).eq(        anti.row(&index) )                 );
        }        
                                                                                                    
        // check that rows and columns are (anti)transposed correctly, even when we pass invalid indices     
        for index in 0 .. 5 {    
            assert!( result_iterators_are_elementwise_equal( 
                matrix.row_result(&index),                 tran.column_result(&index) 
            ));    
            assert!( result_iterators_are_elementwise_equal( 
                matrix.row_result(&index),                 anti.column_reverse_result(&index) 
            ));    

            assert!( result_iterators_are_elementwise_equal( 
                matrix.row_reverse_result(&index),         tran.column_reverse_result(&index) 
            ));  
            assert!( result_iterators_are_elementwise_equal( 
                matrix.row_reverse_result(&index),         anti.column_result(&index) 
            ));                                       
        }

        for index in 0 .. 5 {    
            assert!( result_iterators_are_elementwise_equal( 
                matrix.column_result(&index),                 tran.row_result(&index) 
            ));    
            assert!( result_iterators_are_elementwise_equal( 
                matrix.column_result(&index),                 anti.row_reverse_result(&index) 
            ));    

            assert!( result_iterators_are_elementwise_equal( 
                matrix.column_reverse_result(&index),         tran.row_reverse_result(&index) 
            ));  
            assert!( result_iterators_are_elementwise_equal( 
                matrix.column_reverse_result(&index),         anti.row_result(&index) 
            ));                                       
        }

    }



    



    // Simplify entry-iterators
    // =====================================================================================================================

    // * SEE DOC TESTS    

}




