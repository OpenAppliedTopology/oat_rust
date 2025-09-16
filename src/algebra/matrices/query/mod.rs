//! Provides lookup commands to access sparse matrix entries
//! 
//! This module provides a consistent set of lookup commands that can be used to
//! access rows, columns, and entries in any sparse matrix, regardless of the underlying
//! data structure (CSR, CSC, COO, oracle, etc). These commands are available for any
//! object that implements the [MatrixOracle] trait.
//! 
//! # Content summary
//! - [MatrixOracle] provides lookup commands for rows, columns, and entries
//! - [MatrixAlgebra] provides lookup commands for additional types of data used in algebraic
//!   operations, specifically the coefficient ring of the matrix (encoded as a [ring operator](crate::algebra::rings)),
//!   and the order on row and column indices (encoded as a collection of [order operators](crate::utilities::order)).
//! 
//! # Example
//! 
//! Here's an example of how to use lookup rows, columns, and entries with [MatrixOracle]. Additional examples can be found in
//! the [matrix module](crate::algebra::matrices#build-your-own-matrix).
//! 
//! ```
//! // Import some tools
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec; // a matrix data structure
//! use oat_rust::algebra::matrices::query::MatrixOracle; // trait to retreive entries
//!  
//! // Create a sparse 2x2 upper triangular matrix.
//! let matrix_data =    VecOfVec::new(
//!                         vec![   vec![   (0, 0.),   (1, 1.)  ], 
//!                                 vec![              (1, 2.)  ] 
//!                             ], 
//!                      ).ok().unwrap(); 
//! 
//! // For technical reasons, the matrix oracle traits aren't defined on `VecOfVec` directly.  
//! // Instead, the oracle traits are implemented on references to `VecOfVec`.  References
//! // in Rust are denoted by an `&` symbol.
//! let matrix  =   & matrix_data; 
//! 
//! // Get entries with the `.structural_nonzero_entry` method of the `MatrixOracle` trait
//! assert_eq!( matrix.structural_nonzero_entry( &0, &0 ), Some( 0. ) );
//! assert_eq!( matrix.structural_nonzero_entry( &0, &1 ), Some( 1. ) );
//! assert_eq!( matrix.structural_nonzero_entry( &1, &0 ), None       );
//! assert_eq!( matrix.structural_nonzero_entry( &1, &1 ), Some( 2. ) );
//!  
//! // Look up the nonzero entries in row 0.
//! // To do so, we use the `row` method associated with the `MatrixOracle` trait.
//! let row    =   matrix.row( &0 ); 
//! itertools::assert_equal( row, vec![(0, 0.), (1, 1.)] );
//! 
//! // Look up the nonzero entries in column 1.
//! // To do so, we use the `column` method associated with the `MatrixOracle` trait.
//! let column    =   matrix.column( &1 ); 
//! itertools::assert_equal( column, vec![(0, 1.), (1, 2.)] );
//! ```



pub mod column_helper;




//  =======================================================


use std::fmt::Debug;


use crate::algebra::rings::traits::SemiringOperations;
use crate::algebra::vectors::entries::KeyValGet;
use crate::utilities::order::{JudgeOrder, ReverseOrder};
 // auto-implement a trait on references to objects that implement the trait

 //     =========================================================

        



/// Access sparse matrix entries
/// 
/// 
/// The [MatrixOracle](crate::algebra::matrices::query::MatrixOracle) trait provides a simple, consistent set of commands to retreive information about a sparse matrix.
/// This trait exists because OAT uses many different types of matrices, and we want to use the same set of commands for all of them.
/// 
/// # Examples
/// 
/// - [use this trait](crate::algebra::matrices::query#example)
/// - [make this trait available for a new object](crate::algebra::matrices#build-your-own-matrix)
/// 
/// # Documentation
/// 
/// #### Entries
/// 
/// To get the coefficient in a given row and column, use [`oracle.structural_nonzero_entry( & row, & column )`](crate::algebra::matrices::query::MatrixOracle::structural_nonzero_entry). This 
/// returns `None` if entry `( & row_index, & column_index )` is [structurally zero](#structural-nonzeros-and-invalid-indices), otherwise it returns `Some(coefficient)`.
/// 
/// #### Rows
/// 
/// To access a row, use one of the following 
/// 
/// - [`oracle.row(&index)`](crate::algebra::matrices::query::MatrixOracle::row): Returns a Rust iterator that runs over the entries in the row. Throws an error if the `index` is [invalid](#structural-nonzeros-and-invalid-indices)..
///   See the [documentation](crate::algebra::matrices::query::MatrixOracle::row) for full details.
/// - [`oracle.row_reverse(&index)`](crate::algebra::matrices::query::MatrixOracle::row_reverse): Returns entries in reverse order
/// 
/// #### Columns
/// 
/// To access a column, use one of the following 
/// 
/// - [`oracle.column(&index)`](crate::algebra::matrices::query::MatrixOracle::column): Throws an error if the index is [invalid](#structural-nonzeros-and-invalid-indices).
/// - [`oracle.column_reverse(&index)`](crate::algebra::matrices::query::MatrixOracle::column_reverse): Returns entries in reverse order
/// 
/// #### Error handling for invalid indices
/// 
/// Most look-up methods have corresponding "Result" methods: 
/// [row_result](crate::algebra::matrices::query::MatrixOracle::row_result), 
/// [row_reverse_result](crate::algebra::matrices::query::MatrixOracle::row_reverse_result)
/// [column_result](crate::algebra::matrices::query::MatrixOracle::column_result), 
/// [column_reverse_result](crate::algebra::matrices::query::MatrixOracle::column_reverse_result)
/// etc. The key difference is the way these methods handle [invalid](#structural-nonzeros-and-invalid-indices) inputs; for example, with an [invalid](#structural-nonzeros-and-invalid-indices) row `index`
/// - [matrix_oracl.row(&index)](crate::algebra::matrices::query::MatrixOracle::row) should **panic**
/// - [matrix_oracl.row_result(&index)](crate::algebra::matrices::query::MatrixOracle::row_result) should **not panic** but return a `Result::Err(index)`
/// Returning a `Result::Err(index)` is helpful for debugging and the Rust programming community regards it as a best practice. However, there are several
/// valid reasons why users might avoid using a [Result] method: (1) unfamiliarity -- many users don't know how to use results, (2) readability -- using
/// results can make code a bit harder to read, (3) performance -- because result methods perform extra computation to check validity of inputs, they
/// can consume slightly more computation time. *We have not benchmarked this, but it is possible and even likely that the difference is negligible.
/// in most applied contexts.*
/// 
/// 
/// 
/// #### Structural nonzeros and invalid indices
/// 
/// Here's how we think about some key matrix algebra objects in OAT.
/// - A *dense matrix* is a function `M` from a Cartesian product `I x J` to a set `C`. 
///   - We call `C` the set of *possible coefficients*.
///   - We call `I` the set of *valid row indices*
///   - We call `J` the set of *valid column indices*.
///   - Elements not contained in `I` and `J` are called *invalid*.
/// - A *sparse matrix* is a function `N` from a subset `S` of `I x J` into the set `C`. 
///   - An *invalid entry* of `N` is a pair `(i,j)` such that `i` isn't valid or `j` isn't valid.
///   - A *structural zero entry* of `N` is a pair `(i,j)` such that `i` and `j` are valid, but `S` doesn't contain `(i,j)`
///   - A *structural nonzero entry* of `N` is a tuple `(i,j,c)` such that `(i,j)` is an element of `S` and `N(i,j) =c`
///     - We call `(i,j,c)` a structural nonzero entry *even if `c` is zero*
///     - Every sparse matrix is uniquely determined by the set of it structural nonzero entries.
/// 
/// 
/// Users who implement the [MatrixOracle] trait on a new type of object determine what indices are valid/invalid by implementing
/// the [MatrixOracle::has_row_for_index] and [MatrixOracle::has_column_for_index] methods. They should ensure that the behavior of the
/// lookup commands conforms to the guidelines above.
/// 
/// #### Understanding associated types
/// 
/// There are a lot of different look-up functions, so there are lots of different input and output types.
/// The [MatrixOracle] trait keeps a list of these types, which turns out to useful in some situations.
/// Formally, this bookeeping is handled with [Rust associated types](https://doc.rust-lang.org/rust-by-example/generics/assoc_items/types.html).
/// To help explain what each type means, we've added some "fake code" below.
/// Each line has format 
/// 
/// `object.method_name( InputType) -> OutputType`
/// 
/// Here are some examples (for a complete list, see `Required Methods` and `Provided Methods` in the sidebar):
/// 
/// ```text
/// oracle.structural_nonzero_entry( & RowIndex, & ColumnIndex     )   ->  Option< Coefficient >
/// 
/// oracle.has_row_for_index(     & RowIndex    )     ->  bool
/// oracle.has_column_for_index(  & ColumnIndex )     ->  bool
/// 
/// oracle.row(                     & RowIndex    )     ->  Row
/// oracle.row_result(              & RowIndex    )     ->  Result<Row, RowIndex>   
/// oracle.row_reverse(             & RowIndex    )     ->  RowReverse
/// oracle.row_reverse_restult(     & RowIndex    )     ->  Result<RowReverse, RowIndex>   
/// 
/// oracle.column(                  & ColumnIndex )     ->  Column
/// oracle.column_result(           & ColumnIndex )     ->  Result<Column, ColumnIndex>
/// oracle.column_reverse(          & ColumnIndex )     ->  ColumnReverse
/// oracle.column_reverse_result(   & ColumnIndex )     ->  Result<ColumnReverse, ColumnIndex>
///        
/// row_entry.key()                                     ->  ColumnIndex;
/// row_entry.val()                                     ->  Coefficient;
/// column_entry.key()                                  ->  RowIndex;
/// column_entry.val()                                  ->  Coefficient;
/// ```
/// 
/// 
/// 
/// 
/// # Design notes
/// 
/// There are many different ways to design something like the [MatrixOracle] trait. If you have suggestions for improvement,
/// let us know!  Here are some explanations for some of the design decisions behind the trait
/// 
/// 
/// - Indices and coefficients must implement [Clone], [Debug], and [PartialEq]
///   - So many functions require these conditions that it becomes much easier to
///   enforce them at the trait level, rather than repeat them hundreds of times throughout the
///   library.
/// - Indices require [Eq] but coefficients, row entries, and column entries require only [PartialEq]
///   - The only essential difference between [Eq] and [PartialEq] is that [Eq] requires reflexivity.
///   The float value `NaN` fails this condition, so floats do not implement [Eq]. This one example
///   -- the floating point numbers -- is such an important choice for coefficientst that we 
///   felt we could not require [Eq] for matrix coefficients. However, floats are much less common
///   as indices, and operations like sparse vector addition and multiplication (and by extension,
///   matrix multiplication) could easily yield invaid results if `==` is not reflexive. Thus we impose
///   the stricter [Eq] condition on indices.
/// - Rows and columns must implement [Iterator]
///   - This may appear to introduce some undesirable consequences. For example, it excludes [Vec] as
///   a return type. However, iterating over structural nonzeros is the *only* behavior common to all
///   return types in this library (of which we are currently aware; if this has changed please write
///   to let us know). 
/// 
pub trait MatrixOracle{

    type Coefficient            :   Clone + Debug + PartialEq;    // The type of coefficient stored in each entry of the matrix    
    type RowIndex               :   Clone + Debug + Eq;    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            :   Clone + Debug + Eq;    // The type of column indices    
    type RowEntry               :   Clone + Debug + PartialEq + KeyValGet <  Key = Self::ColumnIndex,  Val = Self::Coefficient   >;  // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            :   Clone + Debug + PartialEq + KeyValGet <  Key = Self::RowIndex,     Val = Self::Coefficient   >;  // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    :   Iterator<   Item    =   Self::RowEntry,     >;  // What you get when you ask for a row.
    type RowReverse             :   Iterator<   Item    =   Self::RowEntry,     >;  // What you get when you ask for a row with the order of entries reversed
    type Column                 :   Iterator<   Item    =   Self::ColumnEntry,  >;  // What you get when you ask for a column
    type ColumnReverse          :   Iterator<   Item    =   Self::ColumnEntry,  >;  // What you get when you ask for a column with the order of entries reversed                             


    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row;
    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse;
    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column;
    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse;
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool;
    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool;

    /// Returns `Some(x)` if there is a structural nonzero at `(row,column)`. 
    /// 
    ///  Returns `None` otherwise.
    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >;

    // The following methods are auto-implemented.
    // ------------------------------------------

    /// Returns `Ok(index)` if `index`` is valid and `Err(index)` if `index` is invalid.
    fn has_row_for_index_result(     &   self, index: & Self::RowIndex   )   -> Result< Self::RowIndex   , Self::RowIndex    >{
        if self.has_row_for_index(index) { Ok(index.clone()) } else { Err(index.clone()) }
    }
    
    /// Returns `Ok(index)` if `index`` is valid and `Err(index)` if `index` is invalid.
    fn has_column_for_index_result(  &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnIndex, Self::ColumnIndex >{
        if self.has_column_for_index(index) { Ok(index.clone()) } else { Err(index.clone()) }
    }

    /// Wraps the output of `structural_nonzero_entry` in a [Result] for error handling.
    /// 
    /// If `self.has_row_for_index(row)` and  `self.has_column_for_index(column)` then returns `Ok(self.structural_nonzero_entry(row,column))`.
    /// Otherwise returns `Err((a,b))`, where `a = self.has_row_for_index_result(row)` and `b = self.has_column_for_index_result(column)`
    fn structural_nonzero_entry_result(               &   self, row:  & Self::RowIndex, column: & Self::ColumnIndex ) 
        ->  
        Result< 
            Option< Self::Coefficient >, 
            ( 
                Result< Self::RowIndex   , Self::RowIndex    >,
                Result< Self::ColumnIndex, Self::ColumnIndex >
            )
        >
    {
        let result_row     =   self.has_row_for_index_result(row);
        let result_col     =   self.has_column_for_index_result(column);
        if result_row.is_err() || result_col.is_err() {
            Err( (result_row, result_col) )
        } else {
            Ok( self.structural_nonzero_entry(row, column) )
        }
    }

    /// Wraps the output of `self.row` in a [Result] for error handling.
    /// 
    /// If `self.has_row_for_index(index)` then returns `Ok(self.row(index))`. Otherwise returns `Err(index)`.
    fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex > 
    { 
        if self.has_row_for_index( index ) { Ok( self.row(index) ) }
        else { Err( index.clone() ) }
    }  

    /// Wraps the output of `self.row_reverse` in a [Result] for error handling.
    /// 
    /// If `self.has_row_for_index(index)` then returns `Ok(self.row(index))`. Otherwise returns `Err(index)`.    
    fn row_reverse_result(         &   self, index: & Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex > 
    { 
        if self.has_row_for_index( index ) { Ok( self.row_reverse(index) ) }
        else { Err( index.clone() ) }
    }

    /// Wraps the output of `self.column` in a [Result] for error handling.
    /// 
    /// If `self.has_column_for_index(index)` then returns `Ok(self.column(index))`. Otherwise returns `Err(index)`.    
    fn column_result(              &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex >
    { 
        if self.has_column_for_index( index ) { Ok( self.column(index) ) }
        else { Err( index.clone() ) }
    }    

    /// Wraps the output of `self.column` in a [Result] for error handling.
    /// 
    /// If `self.has_column_for_index(index)` then returns `Ok(self.column(index))`. Otherwise returns `Err(index)`.      
    fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex >
    { 
        if self.has_column_for_index( index ) { Ok( self.column_reverse(index) ) }
        else { Err( index.clone() ) }
    }        

} 



impl < 'a, T: MatrixOracle >
    
    MatrixOracle for 
    &'a T

{   
    type Coefficient            =   T::Coefficient;       // The type of coefficient stored in each entry of the matrix
    
    type RowIndex               =   T::RowIndex;          // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   T::ColumnIndex;       // The type of column indices

    type RowEntry               =   T::RowEntry   ;          // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   T::ColumnEntry;       // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    
    type Row                    =   T::Row;               // What you get when you ask for a row.
    type RowReverse             =   T::RowReverse;        // What you get when you ask for a row with the order of entries reversed
    type Column                 =   T::Column;            // What you get when you ask for a column   
    type ColumnReverse          =   T::ColumnReverse;     // What you get when you ask for a column with the order of entries reversed                                
    
    fn row(                     &self,  index: &Self::RowIndex    ) -> Self::Row                      { (*self).row( index ) }    
    // fn row_result(                 &   self, index: & Self::RowIndex   )   -> Result< Self::Row, Self::RowIndex >            { (*self).row_result( index ) }        
    fn row_reverse(             &self,  index: &Self::RowIndex    )   -> Self::RowReverse             { (*self).row_reverse( index ) }
    // fn row_reverse_result(         &   self, index: & Self::RowIndex   )   -> Result< Self::RowReverse, Self::RowIndex >     { (*self).row_reverse_result( index ) }    
    
    fn column(                  &self,  index: &Self::ColumnIndex )   -> Self::Column                 { (*self).column( index ) }
    // fn column_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::Column, Self::ColumnIndex >         { (*self).column_result( index ) }
    fn column_reverse(          &self,  index: &Self::ColumnIndex )   -> Self::ColumnReverse          { (*self).column_reverse( index ) }    
    // fn column_reverse_result(      &   self, index: & Self::ColumnIndex)   -> Result< Self::ColumnReverse, Self::ColumnIndex >  { (*self).column_reverse_result( index ) }        

    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool { (*self).has_column_for_index(index) }
    fn has_row_for_index(     &   self, index: &Self::RowIndex   )   -> bool { (*self).has_row_for_index(index) }
    fn structural_nonzero_entry(                   &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
        { (*self).structural_nonzero_entry( row, column ) }
    // fn structural_nonzero_entry_result(               &   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient >
    //     { (*self).structural_nonzero_entry_result( row, column ) }    
} 











//  =======================================================


/// Provides look-up commands for information relevant to matrix algebra
/// 
/// The [MatrixOracle] trait only allows the user to look up entries of a matrix.
/// This isn't enough to perform matrix algebra on it's own.  For example, if we
/// have two matrices with entries of type `usize`, and we want to multiply over
/// the integers modulo 4, then we have to tell the computer that we want to
/// work modulo 4.
/// 
/// The [MatrixAlgebra] trait provides a language to specify this information, via several look-up
/// commands: 
/// 
/// - `matrix.ring_operator()` returns an object called a called an [ring operator](crate::algebra::rings#terminology)
///   that performs algebraic operations on
///   the matrix coefficients
/// - `matrix.order_operator_row_indices()`  returns an
///   object called an [order operator](crate::utilities::order#terminology) that compares the order of row indices.
// ///   For example, if indices are weighted simplices then sometimes you might want to order
// ///   simplices lexicographically, and other times you might want to order them according to weight.
// ///   Whatever your preference, you can define a struct that returns the right value for
// ///   for `matrix.order_operator_row_indices()`. 
///   Similarly, `matrix.order_operator_column_indices()` returns an [order operator](crate::utilities::order#terminology)
///   that compares the order of column indices.
/// - `matrix.order_operator_row_entries()`  returns an [order operator](crate::utilities::order#terminology)
///   for row entries (which are different from row indices). The order on row entries should always agree with column 
///   index order, in the sense that `(i, x) < (j, y)` whenever `i < j`. **NB** because `i` and `j` are
///   column indices, we compare them with the `matrix.order_operator_column_indices()` order operator.
///   You might ask why we need and order operator for entries, if it is supposed to agree with the
///   order operator for column indices.  The answer is that sometimes
///   it can be more efficient to compare two entries directly than to first extract their indices,
///   and then compare those. The `matrix.order_operator_column_entries()` function is similar, but for
///   column indices.
///   
///   **NB** The order that we impose on indices should always match the order in which the
///   [matrix oracle](crate::algebra::matrices#build-your-own-matrix-with-matrix-oracles) returns entries.
///   For example, if `row = matrix.row(&0)` is the iterator returned for row 0 of the matrix, and
///   the first and second entries returned by `row` are `(i,x)` and `(j,y)`, then the [order operator](crate::utilities::order#terminology)
///   for column indices should declare that `i` is strictly less than `j`. The same requirement should
///   be observed for column entries.
/// 
/// 
/// # Where is this used?
/// 
/// OAT uses the [MatrixAlgebra] look-up commands for many of its basic operations: 
/// multiplying matrices, computing dot-products, etc.
/// 
/// # If you have a a  matrix that does note implement [MatrixAlgebra]:
/// 
/// In many cases this is fine; there's a lot you can do without the [MatrixAlgebra] trait.
/// If you do need [MatrixAlgebra] because you want to multiply your matrix with another matrix,
/// or for some other reason, you can wrap your matrix oracle in
/// a [MatrixAlgebraPacket](crate::algebra::matrices::types::packet::MatrixAlgebraPacket).  
pub trait MatrixAlgebra: MatrixOracle {
    type RingOperator:                   Clone + SemiringOperations< Element = Self::Coefficient >;
    type OrderOperatorForRowEntries:     Clone + JudgeOrder< Self::RowEntry >;
    type OrderOperatorForRowIndices:     Clone + JudgeOrder< Self::RowIndex >;
    type OrderOperatorForColumnEntries:  Clone + JudgeOrder< Self::ColumnEntry >;
    type OrderOperatorForColumnIndices:  Clone + JudgeOrder< Self::ColumnIndex >;

    fn ring_operator( &self ) -> Self::RingOperator;
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries;
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices;
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries;
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices;

    fn order_operator_for_row_entries_reverse( self ) -> ReverseOrder< Self::OrderOperatorForRowEntries > 
        where Self: Sized 
    { ReverseOrder::new( self.order_operator_for_row_entries() ) }
    fn order_operator_for_row_indices_reverse( self ) -> ReverseOrder< Self::OrderOperatorForRowIndices > 
        where Self: Sized 
    { ReverseOrder::new( self.order_operator_for_row_indices() ) }
    fn order_operator_for_column_entries_reverse( self ) -> ReverseOrder< Self::OrderOperatorForColumnEntries > 
        where Self: Sized 
    { ReverseOrder::new( self.order_operator_for_column_entries() ) }
    fn order_operator_for_column_indices_reverse( self ) -> ReverseOrder< Self::OrderOperatorForColumnIndices > 
        where Self: Sized 
    { ReverseOrder::new( self.order_operator_for_column_indices() ) } 
}


impl < 'a, T: MatrixAlgebra >
    
    MatrixAlgebra for 
    &'a T
{
    type RingOperator                   =   T::RingOperator;
    type OrderOperatorForRowEntries     =   T::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices     =   T::OrderOperatorForRowIndices;
    type OrderOperatorForColumnEntries  =   T::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices  =   T::OrderOperatorForColumnIndices;

    fn ring_operator( &self ) -> Self::RingOperator { (*self).ring_operator() }
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { (*self).order_operator_for_row_entries() }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { (*self).order_operator_for_row_indices() }    
    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { (*self).order_operator_for_column_entries() }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { (*self).order_operator_for_column_indices() }    
}





