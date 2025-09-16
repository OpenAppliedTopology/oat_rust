//! Memory-efficient scalar matrices (ie diagonal matrices where all entries on the diagonal are equal)

use crate::algebra::matrices::{query::MatrixOracle, operations::MatrixOracleOperations};
use crate::utilities::iterators::general::OncePeekable;
use std::marker::PhantomData;

use std::fmt::Debug;


//  ---------------------------------------------------------------------------
//  SCALAR MATRICES (INDICES CAN BE OF ANY TYPE)
//  ---------------------------------------------------------------------------


//  STRUCT
//  ------

/// Represents a scalar matrix.
///
/// Concretely, this means a diagonal matrix where every diagonal entry equals `alpha`, for some fixed value of `alpha`.
/// 
/// If `m` is a [ScalarMatrix], then the commands `self.row(&index)` and `self.column(&index)`
/// will return an iterator representing the corresponding row or column of the sparse diagonal matrix with
/// `alpha` on the diagaon. Concretely, the return type is `Once< (index, alpha) >`.
/// 
/// This data structure is memory efficient, in the sense that it stores only one copy of `alpha` (rather than one copy of `alpha` for each row of the matrix).
///
/// 
/// # Examples
///
/// ```
/// use oat_rust::algebra::matrices::types::scalar::ScalarMatrix;
/// use oat_rust::algebra::matrices::query::MatrixOracle;   
/// 
/// // create a scalar matrix indexed by `isize` integers
/// let two = ScalarMatrix::new( 2. );
/// assert!(itertools::equal( 
///         two.row( &0 ),   
///         std::iter::once( (0,2.) ),
///     ));  
/// ```
/// 
/// # Design principles for this data structure
/// 
/// It would be possible to write a similar data structure where instead of a `RingOperator`
/// type paratmeter we had a `Coefficient` parameter; the stucture would only need to store
/// a scalar value for the diagonal entries and a different scalar value for the off-diagonal
/// entries. However, because the row and column iterators only return the diagonal elements,
/// it is necessary for the off-diagonal entries to be zero. That means that if we want to 
/// make a "safe" constructor, we want to make it impossible for the user to provide a value
/// for the off-diagonal elements other than zero. Of course "zero" depends on the coefficient
/// ring.
pub struct ScalarMatrix < Index, Scalar >
{
    scalar:         Scalar,
    phantom_key:    PhantomData< Index >
}

impl    < Index, Scalar >
        
        ScalarMatrix 
            < Index, Scalar > 
{
    /// Create new scalar matrix.
    pub fn new( scalar: Scalar ) 
            -> 
            Self  
    {
        ScalarMatrix { 
                scalar,
                phantom_key:        PhantomData,
            }
    }
}


//  ---------------------
//  TRAIT IMPLEMENTATIONS
//  ---------------------


//  MATRIX ORACLE
//  ---------------------------------------------------------------------------

impl    < Index, Scalar >

    MatrixOracle      for 

    ScalarMatrix < Index, Scalar >

    where
        Scalar:     Clone + Debug + PartialEq,
        Index:      Clone + Debug + Eq,     // Clone is needed for (Index,Scalar) to implement KeyValGet
                                    // Eq is needed for the `entry(row,column)` method, to determine if the row and column indices are equal

{
    type Coefficient            =      Scalar; // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =      Index; // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =      Index; // The type of column indices    
    type RowEntry               =      (Index,Scalar);  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =      (Index,Scalar);  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =      OncePeekable< (Index, Scalar) >;  // What you get when you ask for a row.
    type RowReverse             =      OncePeekable< (Index, Scalar) >;  // What you get when you ask for a row with the order of entries reversed
    type Column                 =      OncePeekable< (Index, Scalar) >;  // What you get when you ask for a column
    type ColumnReverse          =      OncePeekable< (Index, Scalar) >;  // What you get when you ask for a column with the order of entries reversed                             

    fn structural_nonzero_entry(                   &   self, row: &Self::RowIndex, column: &Self::ColumnIndex ) -> Option< Self::Coefficient > {
        if row == column { 
            Some( self.scalar.clone() )
        } else { 
            None
        }
    }
    fn has_column_for_index(  &   self, _index: & Self::ColumnIndex)   -> bool 
        { true }
    fn has_row_for_index(     &   self, _index: &Self::RowIndex   )   -> bool 
        { true }
    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row                  
        { OncePeekable::new( ( index.clone(), self.scalar.clone() ) ) } 
    // fn row_result(                 &   self, index: &Self::RowIndex   )   -> Option<Self::Row>
    //     { Some( once( ( index, self.scalar.clone() ) ) ) }  
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse
        { OncePeekable::new( ( index.clone(), self.scalar.clone() ) ) } 
    // fn row_reverse_result(         &   self, index: &Self::RowIndex   )   -> Option<Self::RowReverse>
    //     { Some( once( ( index, self.scalar.clone() ) ) ) } 
    fn column(                  &   self, index: &Self::ColumnIndex)   -> Self::Column
        { OncePeekable::new( ( index.clone(), self.scalar.clone() ) ) }     
    // fn column_result(              &   self, index: &Self::ColumnIndex)   -> Option<Self::Column>
    //     { Some( once( ( index, self.scalar.clone() ) ) ) } 
    fn column_reverse(          &   self, index: &Self::ColumnIndex)   -> Self::ColumnReverse
        { OncePeekable::new( ( index.clone(), self.scalar.clone() ) ) } 
    // fn column_reverse_result(      &   self, index: &Self::ColumnIndex)   -> Option<Self::ColumnReverse>
    //     { Some( once( ( index, self.scalar.clone() ) ) ) } 

} 








//  MATRIX ORACLE OPERATIONS
//  ---------------------------------------------------------------------------

impl    < Index, Scalar >

    MatrixOracleOperations      for 

    ScalarMatrix < Index, Scalar >
{}









//  TESTS
//  =========================================================================================================


//  Doc-test drafts
//  ---------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_test_drafts {
    

    #[test] 
    fn test_scalar_array() {
        use crate::algebra::matrices::types::scalar::ScalarMatrix;
        use crate::algebra::matrices::query::MatrixOracle;   

        // create a scalar matrix indexed by `isize` integers
        let two = ScalarMatrix::<usize, f64>::new( 2. );
        assert!(itertools::equal( 
                two.row( & 0 ),   
                std::iter::once( (0,2.) ),
            ));         
    }    
}










